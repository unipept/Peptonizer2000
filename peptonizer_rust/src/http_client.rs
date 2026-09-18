use std::future::Future;
use std::pin::Pin;

pub type BoxError = Box<dyn std::error::Error + Send + Sync>;
pub type HttpResult<T> = Result<T, BoxError>;
pub type HttpFuture<'a> = Pin<Box<dyn Future<Output = HttpResult<String>> + 'a>>;

pub struct HttpClient {
    #[cfg(not(target_arch = "wasm32"))]
    client: reqwest::Client,
}

impl HttpClient {
    pub fn new() -> Self {
        Self {
            #[cfg(not(target_arch = "wasm32"))]
            client: reqwest::Client::new(),
        }
    }

    pub fn perform_post_request(&self, url: String, payload_json: String) -> HttpFuture<'_> {
        Box::pin(async move {
            #[cfg(target_arch = "wasm32")]
            {
                use js_sys::{Function, Promise, Reflect, global};
                use wasm_bindgen::JsCast;
                use wasm_bindgen_futures::JsFuture;
                use web_sys::{Request, RequestInit, RequestMode, Response};

                let opts = RequestInit::new();
                opts.set_method("POST");
                opts.set_mode(RequestMode::Cors);
                opts.set_body(&payload_json.into());

                let request = Request::new_with_str_and_init(&url, &opts)
                    .map_err(|e| format!("Failed to build request: {:?}", e))?;
                request
                    .headers()
                    .set("Content-Type", "application/json")
                    .map_err(|e| format!("Failed to set request header: {:?}", e))?;

                // Use globalThis.fetch so this works in both Window and Worker contexts.
                let global_this = global();
                let fetch_value = Reflect::get(&global_this, &wasm_bindgen::JsValue::from_str("fetch"))
                    .map_err(|e| format!("Failed to access global fetch function: {:?}", e))?;
                let fetch_fn: Function = fetch_value.dyn_into().map_err(|_| "globalThis.fetch is not callable")?;
                let fetch_result = fetch_fn
                    .call1(&global_this, request.as_ref())
                    .map_err(|e| format!("Failed to invoke fetch: {:?}", e))?;
                let fetch_promise: Promise = fetch_result.dyn_into().map_err(|_| "fetch did not return a Promise")?;

                let response = JsFuture::from(fetch_promise)
                    .await
                    .map_err(|e| format!("Fetch rejected: {:?}", e))?;
                let response: Response = response
                    .dyn_into()
                    .map_err(|e| format!("Failed to decode fetch response: {:?}", e))?;

                if !response.ok() {
                    return Err(format!("Status code {}", response.status()).into());
                }

                let response_text = JsFuture::from(
                    response
                        .text()
                        .map_err(|e| format!("Failed to create response text promise: {:?}", e))?,
                )
                .await
                .map_err(|e| format!("Failed to read response text: {:?}", e))?;
                return response_text
                    .as_string()
                    .ok_or("Response body is not a UTF-8 string".into());
            }

            #[cfg(not(target_arch = "wasm32"))]
            {
                let response = self
                    .client
                    .post(&url)
                    .header("Content-Type", "application/json")
                    .body(payload_json)
                    .send()
                    .await?;

                return Ok(response.text().await?);
            }

            #[allow(unreachable_code)]
            Err("Unsupported target architecture".into())
        })
    }
}

pub fn create_http_client() -> HttpClient {
    HttpClient::new()
}

#[cfg(test)]
mod tests {
    use super::*;
    use serde_json::json;

    #[tokio::test]
    async fn test_invalid_url_returns_error() {
        let client = HttpClient::new();
        let payload = json!({ "message": "hello" });
        let payload = serde_json::to_string(&payload).unwrap();

        let result = client
            .perform_post_request("http://invalid_url".to_string(), payload)
            .await;
        assert!(result.is_err());
    }

    #[tokio::test]
    async fn test_real_http_post() {
        let client = HttpClient::new();
        let payload = json!({ "foo": "bar" });
        let payload = serde_json::to_string(&payload).unwrap();

        let result = client
            .perform_post_request("https://api.unipept.ugent.be".to_string(), payload)
            .await;
        assert!(result.is_ok());
    }
}